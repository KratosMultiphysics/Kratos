###
#
# @copyright (c) 2018 Inria. All rights reserved.
#
###
#
#  @file LibrariesAbsolutePath.cmake
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
#  @date 13-04-2018
#
###

# Transform relative path into absolute path for libraries
# lib_list (input/output): the name of the CMake variable containing libraries, e.g. BLAS_LIBRARIES
# hints_paths (input): additional paths to add when looking for libraries
macro(LIBRARIES_ABSOLUTE_PATH lib_list hints_paths)
  # collect environment paths to dig 
  if(WIN32)
    string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
  elseif(APPLE)
    string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
  else()
    string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
  endif()
  list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
  # copy the lib list 
  set (${lib_list}_COPY "${${lib_list}}")
  # reset the lib list to populate
  set(${lib_list} "")
  foreach(_library ${${lib_list}_COPY})
    if(EXISTS "${_library}")
      # if already an absolute path, nothing special to do
      list(APPEND ${lib_list} ${_library})
    else()
      # replace pattern -lfoo -> foo
      string(REGEX REPLACE "^-l" "" _library "${_library}")
      # remove extensions if exist
      get_filename_component(_ext "${_library}" EXT)
      set(_lib_extensions ".so" ".a" ".dyld" ".dll")
      list(FIND _lib_extensions "${_ext}" _index)
      if (${_index} GREATER -1)
        get_filename_component(_library "${_library}" NAME_WE)
      endif()
      # try to find the lib
      find_library(_library_path NAMES ${_library} HINTS ${hints_paths} ${_lib_env})
      if (_library_path)
          list(APPEND ${lib_list} ${_library_path})
      else()
          message(FATAL_ERROR "Dependency of ${lib_list} '${_library}' NOT FOUND")
      endif()
      unset(_library_path CACHE)
    endif()
  endforeach()
endmacro()

##
## @end file LibrariesAbsolutePath.cmake
##
