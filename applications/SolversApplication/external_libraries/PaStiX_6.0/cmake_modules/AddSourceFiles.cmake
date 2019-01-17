#
# See https://cmake.org/pipermail/cmake/2010-March/035992.html
#
function(add_source_files)
  cmake_parse_arguments(add_source_files "" "DIRECTORY" "" ${ARGN})

  get_property(is_defined GLOBAL PROPERTY PASTIX_LIB_SRCS2 DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY PASTIX_LIB_SRCS2
      BRIEF_DOCS "List of source files included in the PaStiX library"
      FULL_DOCS "List of source files included in the PaStiX library")
  endif()

  if(NOT ${add_source_files_DIRECTORY} STREQUAL "")
    set(_dir "${add_source_files_DIRECTORY}/")
  else()
    set(_dir "")
  endif()

  set(ARGN ${add_source_files_UNPARSED_ARGUMENTS})
  # make absolute paths
  set(SRCS)
  foreach(s IN LISTS ARGN)
    set(s "${_dir}${s}")
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${s}")
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY PASTIX_LIB_SRCS2 "${SRCS}")
endfunction(add_source_files)

#
# See https://cmake.org/pipermail/cmake/2010-March/035992.html
#
function(add_documented_files)
  cmake_parse_arguments(add_documented_files "" "DIRECTORY" "" ${ARGN})

  get_property(is_defined GLOBAL PROPERTY PASTIX_DOX_SRCS DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY PASTIX_DOX_SRCS
      BRIEF_DOCS "List of documented source files"
      FULL_DOCS "List of source files to be included into the in-code documentation")
  endif()

  if(NOT ${add_documented_files_DIRECTORY} STREQUAL "")
    set(_dir "${add_documented_files_DIRECTORY}/")
  else()
    set(_dir "")
  endif()

  set(ARGN ${add_documented_files_UNPARSED_ARGUMENTS})
  # make absolute paths
  set(SRCS)
  foreach(s IN LISTS ARGN)
    set(s "${_dir}${s}")
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${s}")
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY PASTIX_DOX_SRCS "${SRCS}")
endfunction(add_documented_files)
