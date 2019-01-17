###
#
# @copyright (c) 2018 Inria. All rights reserved.
#
###
#
#  @file FindPkgconfigLibrariesAbsolutePath.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 1.0.0
#  @author Florent Pruvost
#  @date 06-04-2018
#
###

# Transform relative path into absolute path for libraries found with the
# pkg_search_module cmake macro
# _prefix: the name of the CMake variable used when pkg_search_module was called
# e.g. for pkg_search_module(BLAS blas) _prefix would be BLAS
macro(FIND_PKGCONFIG_LIBRARIES_ABSOLUTE_PATH _prefix)
  if(WIN32)
    string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
  elseif(APPLE)
    string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
  else()
    string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
  endif()
  list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
  # non static case
  set(${_prefix}_LIBRARIES_COPY "${${_prefix}_LIBRARIES}")
  set(${_prefix}_LIBRARIES "")
  foreach(_library ${${_prefix}_LIBRARIES_COPY})
    if(EXISTS "${_library}")
      list(APPEND ${_prefix}_LIBRARIES ${_library})
    else()
      get_filename_component(_ext "${_library}" EXT)
      set(_lib_extensions ".so" ".a" ".dyld" ".dll")
      list(FIND _lib_extensions "${_ext}" _index)
      if (${_index} GREATER -1)
        get_filename_component(_library "${_library}" NAME_WE)
      endif()
      find_library(_library_path NAMES ${_library}
          HINTS ${${_prefix}_LIBDIR} ${${_prefix}_LIBRARY_DIRS} ${_lib_env})
      if (_library_path)
          list(APPEND ${_prefix}_LIBRARIES ${_library_path})
      else()
          message(FATAL_ERROR "Dependency of ${_prefix} '${_library}' NOT FOUND")
      endif()
      unset(_library_path CACHE)
    endif()
  endforeach()
  set (${_prefix}_LIBRARIES "${${_prefix}_LIBRARIES}" CACHE INTERNAL "" FORCE)
  ## static case
  #set(${_prefix}_STATIC_LIBRARIES_COPY "${${_prefix}_STATIC_LIBRARIES}")
  #set(${_prefix}_STATIC_LIBRARIES "")
  #foreach(_library ${${_prefix}_STATIC_LIBRARIES_COPY})
  #  if(EXISTS "${_library}")
  #    list(APPEND ${_prefix}_STATIC_LIBRARIES ${_library})
  #  else()
  #    get_filename_component(_ext "${_library}" EXT)
  #    set(_lib_extensions ".so" ".a" ".dyld" ".dll")
  #    list(FIND _lib_extensions "${_ext}" _index)
  #    if (${_index} GREATER -1)
  #      get_filename_component(_library "${_library}" NAME_WE)
  #    endif()
  #    # try static first
  #    set (default_find_library_suffixes ${CMAKE_FIND_LIBRARY_SUFFIXES})
  #    set (CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
  #    find_library(_library_path NAMES ${_library}
  #        HINTS ${${_prefix}_STATIC_LIBDIR} ${${_prefix}_STATIC_LIBRARY_DIRS} ${_lib_env})
  #    set (CMAKE_FIND_LIBRARY_SUFFIXES ${default_find_library_suffixes})
  #    # if not found try dynamic
  #    if (NOT _library_path)
  #      find_library(_library_path NAMES ${_library}
  #          HINTS ${${_prefix}_STATIC_LIBDIR} ${${_prefix}_STATIC_LIBRARY_DIRS} ${_lib_env})
  #    endif()
  #    if (_library_path)
  #        list(APPEND ${_prefix}_STATIC_LIBRARIES ${_library_path})
  #    else()
  #        message(FATAL_ERROR "Dependency of ${_prefix} '${_library}' NOT FOUND")
  #    endif()
  #    unset(_library_path CACHE)
  #  endif()
  #endforeach()
  #set (${_prefix}_STATIC_LIBRARIES "${${_prefix}_STATIC_LIBRARIES}" CACHE INTERNAL "" FORCE)
endmacro()

##
## @end file FindPkgconfigLibrariesAbsolutePath.cmake
##
