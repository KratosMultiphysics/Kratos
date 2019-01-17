###
#
# @copyright (c) 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                          Univ. Bordeaux. All rights reserved.
#
###
#
#  @file GenPkgConfig.cmake
#
#  @project MORSE
#
#  @version 0.9.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2018-06-21
#
#.rst:
# GenPkgConfig
# ------------
#
#   The GenPkgConfig module intend to provide functions to generate a
#   pkg-config file with the right format, and a environment file that
#   the user can source to setup its environment with the library
#   installed.
#
# generate_pkgconfig_files(
#   file1 file2 ...
#   [PROJECTNAME name]
#   [LIBS lib1 lib2 ...]
#   [LIBS_PRIVATE lib1 lib2 ...]
#   [REQUIRED pkg1 pkg2 ...]
#   [REQUIRED_PRIVATE pkg1 pkg2 ...]
# )
# Where:
#   - file1, file2, .. are the different input files for the
#     pkg-config configuration files.
#   - PROJECTNAME defines the prefix of the variables to use for the
#     default following values ${PROJECTNAME}_PKGCONFIG_XXX with XXX
#     being part of (LIBS, LIBS_PRIVATE, REQUIRED and
#     REQUIRED_PRIVATE)
#   - LIBS defines the libs of the current project for dynamic
#     linking.
#   - LIBS_PRIVATE defines the libs of the current project and its
#     required dependencies not handled by pkg-config for static
#     linking.
#   - LIBS defines the additional pkg-config packages needed by the
#     project for dynamic linking.
#   - LIBS_PRIVATE defines the additional pkg-config packages needed
#     by the project for static linking.
#
# The generated files are then installed in the subdirectory
# lib/pkgconfig of the prefix directory.
#
# generate_env_files(
#   [PROJECTNAME name]
# )
# Where:
#   - PROJECTNAME defines the project name of the library if different
#     from CMAKE_PROJECT_NAME.
#
# The generated file is installed in the subdirectory bin of the
# prefix directory.
#
###

###
#
# gpc_convert_incstyle_to_pkgconfig(): convert a libraries list to
# follow the pkg-config style used in gpc_clean_lib_list
#
###
macro(gpc_convert_incstyle_to_pkgconfig _inclist)
  set(${_inclist}_CPY "${${_inclist}}")
  set(${_inclist} "")
  foreach(_dep ${${_inclist}_CPY})
    if (${_dep} MATCHES "^-D")
      list(APPEND ${_inclist} ${_dep})
    else()
      list(APPEND ${_inclist} "-I${_dep}")
    endif()
  endforeach()
endmacro(gpc_convert_incstyle_to_pkgconfig)

###
#
# gpc_convert_libstyle_to_pkgconfig(): convert a libraries list to
# follow the pkg-config style used in gpc_clean_lib_list
#
###
macro(gpc_convert_libstyle_to_pkgconfig _liblist)

  # Start by removing duplicates
  list(REVERSE ${_liblist})
  list(REMOVE_DUPLICATES ${_liblist})
  list(REVERSE ${_liblist})

  # Convert to pkg-config file format
  set(${_liblist}_CPY "${${_liblist}}")
  set(${_liblist} "")
  foreach(_dep ${${_liblist}_CPY})
    if (${_dep} MATCHES "^/")
      get_filename_component(dep_libname ${_dep} NAME)
      get_filename_component(dep_libdir  ${_dep} PATH)
      string(REPLACE ".so"    "" dep_libname "${dep_libname}")
      string(REPLACE ".a"     "" dep_libname "${dep_libname}")
      string(REPLACE ".dylib" "" dep_libname "${dep_libname}")
      string(REPLACE ".dll"   "" dep_libname "${dep_libname}")
      string(REPLACE "lib"    "" dep_libname "${dep_libname}")
      list(APPEND ${_liblist} -L${dep_libdir} -l${dep_libname})
    elseif(NOT ${_dep} MATCHES "^-")
      list(APPEND ${_liblist} "-l${_dep}")
    else()
      list(APPEND ${_liblist} ${_dep})
    endif()
  endforeach()
endmacro(gpc_convert_libstyle_to_pkgconfig)

###
#
# gpc_clean_lib_list(): clean libraries lists to follow the pkg-config
# style used in GENERATE_PKGCONFIG_FILE
#
###
macro(gpc_clean_lib_list _package)
  if ( ${_package}_PKGCONFIG_INCS )
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_INCS)
    gpc_convert_incstyle_to_pkgconfig(${_package}_PKGCONFIG_INCS)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_INCS "${${_package}_PKGCONFIG_INCS}")
  endif()
  if ( ${_package}_PKGCONFIG_LIBS )
    gpc_convert_libstyle_to_pkgconfig(${_package}_PKGCONFIG_LIBS)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS "${${_package}_PKGCONFIG_LIBS}")
  endif()
  if ( ${_package}_PKGCONFIG_LIBS_PRIVATE )
    gpc_convert_libstyle_to_pkgconfig(${_package}_PKGCONFIG_LIBS_PRIVATE)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS_PRIVATE "${${_package}_PKGCONFIG_LIBS_PRIVATE}")
  endif()
  if ( ${_package}_PKGCONFIG_REQUIRED )
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED "${${_package}_PKGCONFIG_REQUIRED}")
  endif()
  if ( ${_package}_PKGCONFIG_REQUIRED_PRIVATE )
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED_PRIVATE)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED_PRIVATE "${${_package}_PKGCONFIG_REQUIRED_PRIVATE}")
  endif()
endmacro(gpc_clean_lib_list)

###
#
# generate_pkgconfig_file: generate pkg-config file
#
###
macro(generate_pkgconfig_files)
  set(_options )
  set(_oneValueArgs PROJECTNAME)
  set(_multiValueArgs LIBS LIBS_PRIVATE REQUIRED REQUIRED_PRIVATE)
  cmake_parse_arguments(generate_pkgconfig_file
    "${_options}" "${_oneValueArgs}"
    "${_multiValueArgs}" ${ARGN} )

  if ( NOT DEFINED PROJECTNAME )
    set( PROJECTNAME ${CMAKE_PROJECT_NAME} )
  endif()

  set(ARGN ${generate_pkgconfig_file_UNPARSED_ARGUMENTS})

  # The link flags specific to this package and any required libraries
  # that don't support PkgConfig
  list(APPEND ${PROJECTNAME}_PKGCONFIG_LIBS ${LIBS})

  # The link flags for private libraries required by this package but not
  # exposed to applications
  list(APPEND ${PROJECTNAME}_PKGCONFIG_LIBS_PRIVATE ${LIBS_PRIVATE})

  # A list of packages required by this package
  list(APPEND ${PROJECTNAME}_PKGCONFIG_REQUIRED ${REQUIRED})

  # A list of private packages required by this package but not exposed to
  # applications
  list(APPEND ${PROJECTNAME}_PKGCONFIG_REQUIRED_PRIVATE ${REQUIRED_PRIVATE})

  # Define required package
  # -----------------------
  gpc_clean_lib_list(${PROJECTNAME})

  foreach(f IN LISTS ARGN)
    get_filename_component(fname "${f}" NAME_WE)

    # Create .pc files
    # ----------------
    configure_file(
      ${f}
      "${CMAKE_BINARY_DIR}/lib/pkgconfig/${fname}.pc" @ONLY)

    # installation
    # ------------
    install(FILES
      "${CMAKE_BINARY_DIR}/lib/pkgconfig/${fname}.pc"
      DESTINATION lib/pkgconfig )

  endforeach()

endmacro(generate_pkgconfig_files)

###
#
# generate_env_file: generate pkf-config files
#
###
macro(generate_env_file)

  set(_options )
  set(_oneValueArgs PROJECTNAME)
  set(_multiValueArgs )
  cmake_parse_arguments(generate_env_file
    "${_options}" "${_oneValueArgs}"
    "${_multiValueArgs}" ${ARGN} )

  if ( NOT DEFINED PROJECTNAME )
    set( PROJECTNAME ${CMAKE_PROJECT_NAME} )
  endif()

  string(TOLOWER ${PROJECTNAME} LONAME)
  string(TOUPPER ${PROJECTNAME} UPNAME)

  # Create .sh file
  # ---------------
  configure_file(
    "${MORSE_CMAKE_MODULE_PATH}/env.sh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/bin/${LONAME}_env.sh" @ONLY)

  # installation
  # ------------
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/bin/${LONAME}_env.sh"
    DESTINATION bin)

endmacro(generate_env_file)

##
## @end file GenPkgConfig.cmake
##
