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

# Create compilation modes:
#     - Maintener: adds compiler warnings to Debug mode
#     - RelWithAssert: RelWithDebInfo without NDEBUG precompiler flag

# Add Maintener mode
if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set ( CLANG_FLAGS
    "-O0 -g -Weverything -Wno-sign-conversion -Wno-char-subscripts -Wno-padded")
  set(CMAKE_CXX_FLAGS_MAINTAINER "${CLANG_FLAGS}"
    CACHE STRING
    "Flags used by the CXX compiler during Maintainer builds."
    FORCE)
  set(CMAKE_C_FLAGS_MAINTAINER "${CLANG_FLAGS}"
    CACHE STRING
    "Flags used by the C compiler during Maintainer builds."
    FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_MAINTAINER ${CMAKE_EXE_LINKER_FLAGS_DEBUG}
    CACHE STRING
    "Flags used for linking binaries during Maintainer builds."
    FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS_MAINTAINER ${CMAKE_SHARED_LINKER_FLAGS_DEBUG}
    CACHE STRING
    "Flags used by the shared libraries linker during Maintainer builds."
    FORCE)
  set(CMAKE_STATIC_LINKER_FLAGS_MAINTAINER ${CMAKE_STATIC_LINKER_FLAGS_DEBUG}
    CACHE STRING
    "Flags used by the static libraries linker during Maintainer builds."
    FORCE)

elseif (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_MAINTAINER "-O0 -g -Wall" CACHE STRING
    "Flags used by the CXX compiler during Maintainer builds."
    FORCE)
  set(CMAKE_C_FLAGS_MAINTAINER "-O0 -g -Wall" CACHE STRING
    "Flags used by the C compiler during Maintainer builds."
    FORCE)

  set ( LD_LINKER_FLAGS "-Wl,--warn-unresolved-symbols,--warn-once" )

  set(CMAKE_EXE_LINKER_FLAGS_MAINTAINER
    ${LD_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}
    CACHE STRING
    "Flags used for linking binaries during Maintainer builds."
    FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
    ${LD_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_DEBUG}
    CACHE STRING
    "Flags used by the shared libraries linker during Maintainer builds."
    FORCE)
  # Static lib linking uses ar and not ld: -Wl is not supported
  set(CMAKE_STATIC_LINKER_FLAGS_MAINTAINER
    ${CMAKE_STATIC_LINKER_FLAGS_DEBUG}
    CACHE STRING
    "Flags used by the static libraries linker during Maintainer builds."
    FORCE)

else ()
  # Not implemented: use Debug flags
  set(CMAKE_CXX_FLAGS_MAINTAINER "${CMAKE_CXX_FLAGS_DEBUG}"
    CACHE STRING
    "Flags used by the CXX compiler during Maintainer builds."
    FORCE)
  set(CMAKE_C_FLAGS_MAINTAINER "${CMAKE_C_FLAGS_DEBUG}"
    CACHE STRING
    "Flags used by the C compiler during Maintainer builds."
    FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_MAINTAINER "${CMAKE_EXE_LINKER_FLAGS_DEBUG}"
    CACHE STRING
    "Flags used for linking binaries during Maintainer builds."
    FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS_MAINTAINER "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}"
    CACHE STRING
    "Flags used by the shared libraries linker during Maintainer builds."
    FORCE)
  set(CMAKE_STATIC_LINKER_FLAGS_MAINTAINER "${CMAKE_STATIC_LINKER_FLAGS_DEBUG}"
    CACHE STRING
    "Flags used by the static libraries linker during Maintainer builds."
    FORCE)

endif()

mark_as_advanced(
  CMAKE_CXX_FLAGS_MAINTAINER
  CMAKE_C_FLAGS_MAINTAINER
  CMAKE_EXE_LINKER_FLAGS_MAINTAINER
  CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
  CMAKE_STATIC_LINKER_FLAGS_MAINTAINER
  )

# Add RelWithAssert mode
STRING(REGEX REPLACE ".DNDEBUG" " "
  RELWITHASSERT_C_FLAGS "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )

STRING(REGEX REPLACE ".DNDEBUG" " "
  RELWITHASSERT_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" )

set(CMAKE_CXX_FLAGS_RELWITHASSERT "${RELWITHASSERT_CXX_FLAGS}"
  CACHE STRING
  "Flags used by the CXX compiler during RelWithAssert builds."
  FORCE)

set(CMAKE_C_FLAGS_RELWITHASSERT "${RELWITHASSERT_C_FLAGS}"
  CACHE STRING
  "Flags used by the C compiler during RelWithAssert builds."
  FORCE)

set(CMAKE_EXE_LINKER_FLAGS_RELWITHASSERT "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}"
  CACHE STRING
  "Flags used for linking binaries during RelWithAssert builds."
  FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_RELWITHASSERT "${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEINFO}"
  CACHE STRING
  "Flags used by the shared libraries linker during maintainer builds."
  FORCE)
set(CMAKE_STATIC_LINKER_FLAGS_RELWITHASSERT "${CMAKE_STATIC_LINKER_FLAGS_RELWITHDEINFO}"
  CACHE STRING
  "Flags used by the static libraries linker during maintainer builds."
  FORCE)

mark_as_advanced(
  CMAKE_CXX_FLAGS_RELWITHASSERT
  CMAKE_C_FLAGS_RELWITHASSERT
  CMAKE_EXE_LINKER_FLAGS_RELWITHASSERT
  CMAKE_SHARED_LINKER_FLAGS_RELWITHASSERT
  CMAKE_STATIC_LINKER_FLAGS_RELWITHASSERT
  )
