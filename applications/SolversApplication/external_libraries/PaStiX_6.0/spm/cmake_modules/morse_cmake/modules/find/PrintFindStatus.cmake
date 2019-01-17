###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Some macros to print status when search for headers and libs
# Main parameters of macros
#  _libname: name of the lib you seek, foo for example
#  _header_to_find: name of the header you seek, foo.h for example
#  _lib_to_find: name of the library you seek, libfoo for example
#  _pc_to_find: name of the pkg-config file zyou seek, foo.pc for example


#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.

#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)


# Set some colors
#if(NOT WIN32)
#    string(ASCII 27 Esc)
#    set(ColourReset "${Esc}[m")
#    set(ColourBold  "${Esc}[1m")
#    set(Red         "${Esc}[31m")
#    set(Green       "${Esc}[32m")
#    set(Yellow      "${Esc}[33m")
#    set(Blue        "${Esc}[34m")
#    set(Magenta     "${Esc}[35m")
#    set(Cyan        "${Esc}[36m")
#    set(White       "${Esc}[37m")
#    set(BoldRed     "${Esc}[1;31m")
#    set(BoldGreen   "${Esc}[1;32m")
#    set(BoldYellow  "${Esc}[1;33m")
#    set(BoldBlue    "${Esc}[1;34m")
#    set(BoldMagenta "${Esc}[1;35m")
#    set(BoldCyan    "${Esc}[1;36m")
#    set(BoldWhite   "${Esc}[1;37m")
#endif()


# This macro informs why the _header_to_find file has not been found
macro(Print_Find_Header_Status _libname _header_to_find)

  # save _libname upper and lower case
  string(TOUPPER ${_libname} LIBNAME)
  string(TOLOWER ${_libname} libname)

  # print status
  #message(" ")
  if(${LIBNAME}_INCDIR)
    message("${Blue}${LIBNAME}_INCDIR is defined but ${_header_to_find}"
      "has not been found in ${${LIBNAME}_INCDIR}${ColourReset}")
  else()
    if(${LIBNAME}_DIR)
      message("${Blue}${LIBNAME}_DIR is defined but"
	"${_header_to_find} has not been found in"
	"${${LIBNAME}_DIR}/include${ColourReset}")
    else()
      message("${Blue}${_header_to_find} not found."
	"Nor ${LIBNAME}_DIR neither ${LIBNAME}_INCDIR"
	"are defined so that we looked for ${_header_to_find} in"
	"system paths (INCLUDE, CPATH, C_INCLUDE_PATH,"
	"INCLUDE_PATH, CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES)${ColourReset}")
      if(_inc_env)
	message("${Blue}${_header_to_find} has not been found in"
	  "${_inc_env}${ColourReset}")
      endif()
    endif()
  endif()
  message("${BoldBlue}Please indicate where to find ${_header_to_find}. You have three options:\n"
    "- Option 1: Provide the root directory of the library with cmake option: -D${LIBNAME}_DIR=your/path/to/${libname}/\n"
    "- Option 2: Provide the directory where to find the headers with cmake option: -D${LIBNAME}_INCDIR=your/path/to/${libname}/include/\n"
    "- Option 3: Update your environment variable (INCLUDE or CPATH)\n"
    "- Option 4: If your library provides a PkgConfig file, make sure pkg-config finds your library${ColourReset}")
  #message(" ")

endmacro()

# This macro informs why the _lib_to_find file has not been found
macro(Print_Find_Library_Status _libname _lib_to_find)

  # save _libname upper/lower case
  string(TOUPPER ${_libname} LIBNAME)
  string(TOLOWER ${_libname} libname)

  # print status
  #message(" ")
  if(${LIBNAME}_LIBDIR)
    message("${Yellow}${LIBNAME}_LIBDIR is defined but ${_lib_to_find}"
      "has not been found in ${${LIBNAME}_LIBDIR}${ColourReset}")
  else()
    if(${LIBNAME}_DIR)
      message("${Yellow}${LIBNAME}_DIR is defined but ${_lib_to_find}"
	"has not been found in ${${LIBNAME}_DIR}/lib(or /lib32 or"
	"/lib64)${ColourReset}")
    else()
      message("${Yellow}${_lib_to_find} not found."
	"Nor ${LIBNAME}_DIR neither ${LIBNAME}_LIBDIR"
	"are defined so that we looked for ${_lib_to_find} in"
	"system paths (Linux: LD_LIBRARY_PATH, Windows: LIB,"
	"Mac: DYLD_LIBRARY_PATH,"
	"CMAKE_C_IMPLICIT_LINK_DIRECTORIES)${ColourReset}")
      if(_lib_env)
	message("${Yellow}${_lib_to_find} has not been found in"
	  "${_lib_env}${ColourReset}")
      endif()
    endif()
  endif()
  message("${BoldYellow}Please indicate where to find ${_lib_to_find}. You have three options:\n"
    "- Option 1: Provide the root directory of the library with cmake option: -D${LIBNAME}_DIR=your/path/to/${libname}/\n"
    "- Option 2: Provide the directory where to find the library with cmake option: -D${LIBNAME}_LIBDIR=your/path/to/${libname}/lib/\n"
    "- Option 3: Update your environment variable (Linux: LD_LIBRARY_PATH, Windows: LIB, Mac: DYLD_LIBRARY_PATH)\n"
    "- Option 4: If your library provides a PkgConfig file, make sure pkg-config finds your library${ColourReset}")

endmacro()

# This macro informs why the _lib_to_find file has not been found
macro(Print_Find_Library_Blas_Status _libname _lib_to_find)

  # save _libname upper/lower case
  string(TOUPPER ${_libname} LIBNAME)
  string(TOLOWER ${_libname} libname)

  # print status
  #message(" ")
  if(${LIBNAME}_LIBDIR)
    message("${Yellow}${LIBNAME}_LIBDIR is defined but ${_lib_to_find}"
      "has not been found in ${ARGN}${ColourReset}")
  else()
    if(${LIBNAME}_DIR)
      message("${Yellow}${LIBNAME}_DIR is defined but ${_lib_to_find}"
	"has not been found in ${ARGN}${ColourReset}")
    else()
      message("${Yellow}${_lib_to_find} not found."
	"Nor ${LIBNAME}_DIR neither ${LIBNAME}_LIBDIR"
	"are defined so that we look for ${_lib_to_find} in"
	"system paths (Linux: LD_LIBRARY_PATH, Windows: LIB,"
	"Mac: DYLD_LIBRARY_PATH,"
	"CMAKE_C_IMPLICIT_LINK_DIRECTORIES)${ColourReset}")
      if(_lib_env)
	message("${Yellow}${_lib_to_find} has not been found in"
	  "${_lib_env}${ColourReset}")
      endif()
    endif()
  endif()
  message("${BoldYellow}Please indicate where to find ${_lib_to_find}. You have three options:\n"
    "- Option 1: Provide the root directory of the library with cmake option: -D${LIBNAME}_DIR=your/path/to/${libname}/\n"
    "- Option 2: Provide the directory where to find the library with cmake option: -D${LIBNAME}_LIBDIR=your/path/to/${libname}/lib/\n"
    "- Option 3: Update your environment variable (Linux: LD_LIBRARY_PATH, Windows: LIB, Mac: DYLD_LIBRARY_PATH)\n"
    "- Option 4: If your library provides a PkgConfig file, make sure pkg-config finds your library${ColourReset}")

endmacro()

# This macro informs why the _lib_to_find file has not been found
macro(Print_Find_Library_Blas_CheckFunc_Status _name)

  # save _libname upper/lower case
  string(TOUPPER ${_name} FUNCNAME)
  string(TOLOWER ${_name} funcname)

  # print status
  #message(" ")
  message("${Red}Libs have been found but check of symbol ${_name} failed "
    "with following libraries ${ARGN}${ColourReset}")
  message("${BoldRed}Please open your error file CMakeFiles/CMakeError.log"
    "to figure out why it fails${ColourReset}")
  #message(" ")

endmacro()

# This macro informs that _pc_to_find file has not been found in the list
# path you give as last argument (read in ${ARGN})
# ex: Print_Find_Pkgconfig_Status(foo foo.pc ${PATHLIST}
macro(Print_Find_Pkgconfig_Status _libname _pc_to_find)

  # save _libname lower case
  string(TOLOWER ${_libname} libname)

  # print status
  #message(" ")
  message("${Magenta}${_pc_to_find} has not been found in"
    "${ARGN}${ColourReset}")
  message("${BoldMagenta}If you really want to use the pkg-config file of"
    "${libname}, please update your PKG_CONFIG_PATH with the path"
    "where ${_pc_to_find} states${ColourReset}")
  #message(" ")

endmacro()
